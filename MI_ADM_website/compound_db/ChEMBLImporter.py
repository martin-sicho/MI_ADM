from chembl_webresource_client import CompoundResource
from chembl_webresource_client import TargetResource
from compound_db import models
from django.conf import settings
from django.db import transaction, IntegrityError
from django.utils import timezone
from rdkit import Chem

from MI_ADM_website.compound_db.utils import is_number


class ChEMBLError(Exception):

    def __init__(self, compound_id, error_msg):
        super(ChEMBLError, self).__init__('Failed to fetch ChEMBL data for {0} ({1})'.format(compound_id, error_msg))

class ChEMBLImporter:

    CHEMBL_TARGETS_RESOURCE = TargetResource()
    CHEMBL_COMPOUNDS_RESOURCE = CompoundResource()
    # TODO: test if connection OK before attempting retrieval

    VALUE_OPERATORS_MAP = {
        '=' : lambda thrs, val : val == thrs
        , '>=' : lambda thrs, val : val >= thrs
        , '>' : lambda thrs, val : val > thrs
        , '<=' : lambda thrs, val : val <= thrs
        , '<' : lambda thrs, val : val < thrs
    }

    def __init__(self, target_id, description, filters=None):
        self.skipped_compounds = [] # TODO: log this in the database for future review
        self.filtered_compounds = []
        self.target = self.CHEMBL_TARGETS_RESOURCE.get(target_id)
        self.description = description
        self.is_data_saved = False
        self.exceptions = []
        self.fatal_exception = None
        self.filters = filters

    def _check_response(self, compound_data, activity_info):
        if type(compound_data) != dict:
            if type(compound_data) == int:
                self.exceptions.append(ChEMBLError(activity_info['ingredient_cmpd_chemblid'], compound_data))
            else:
                self.exceptions.append(ChEMBLError(activity_info['ingredient_cmpd_chemblid'], 'Did not get a dictionary'))
            return False
        if 'smiles' not in compound_data:
            self.exceptions.append(ChEMBLError(activity_info['ingredient_cmpd_chemblid'], 'No SMILES found'))
            return False
        if not is_number(activity_info['value']):
            self.exceptions.append(ChEMBLError(activity_info['ingredient_cmpd_chemblid'], 'Activity value cannot be converted to a single number'))
            return False
        return True

    def _check_filter(self, key, activity_info):
        if self.filters[key] \
                and 'All' not in self.filters[key] \
                and activity_info[key] not in self.filters[key]:
            return False
        else:
            return True

    def _apply_filters(self, activity_info):
        status = True
        special = {'value', 'activity_comment'}
        for key in self.filters:
            if key == 'value' and self.filters[key][0]:
                thrs = self.filters[key][0]
                oper = self.filters[key][1]
                activity_value = float(activity_info[key]) if is_number(activity_info[key]) else None
                if activity_value:
                    status = self.VALUE_OPERATORS_MAP[oper](thrs, activity_value)
            elif key == 'activity_comment':
                self.filters[key].discard('')
                comment = activity_info[key].strip()
                if self.filters[key] and comment not in self.filters[key]:
                    status = False
            elif key not in special:
                status = self._check_filter(key, activity_info)
            if not status:
                break
        return status

    def save_data(self):
        target_info = self.target
        try:
            with transaction.atomic():
                target_data = models.ChEMBLTargetData(
                    uniprot_accession=target_info['proteinAccession']
                    , chembl_id=target_info['chemblId']
                    , organism=target_info['organism']
                    , preffered_name=target_info['preferredName']
                    , description=self.description
                )
                target_data.save()

                activities = self.CHEMBL_TARGETS_RESOURCE.bioactivities(target_info['chemblId'])
                for idx,activity_info in enumerate(activities):
                    # print('Processing {0}/{1}...'.format(idx, len(activities)))

                    if not self._apply_filters(activity_info):
                        self.filtered_compounds.append(activity_info)
                        continue

                    compound_data = self.CHEMBL_COMPOUNDS_RESOURCE.get(activity_info['ingredient_cmpd_chemblid'])
                    if not self._check_response(compound_data, activity_info):
                        self.skipped_compounds.append(compound_data)
                        continue

                    mol = Chem.MolFromSmiles(compound_data['smiles'])
                    if mol:
                        compound = models.Compound(mol, 'Imported from ChEMBL on {0}'.format(timezone.now()))
                        try:
                            with transaction.atomic():
                                compound.save()
                        except models.Compound.MoleculeAlreadyExists as exp:
                            self.exceptions.append(exp)
                            compound = models.Compound.objects.get(inchi_key=compound.inchi_key)
                        except IntegrityError as ierr:
                            self.exceptions.append(ierr)
                            self.skipped_compounds.append(compound_data) # FIXME: this results in a loss of activity data if an enantiomer is already in the database
                            continue

                        bioassay_data = models.ChEMBLBioassayData(
                            compound=compound
                            , target_data=target_data
                            , assay_id=activity_info['assay_chemblid']
                            , ingredient_cmpd_id=activity_info['ingredient_cmpd_chemblid']
                            , units=activity_info['units']
                            , bioactivity_type=activity_info['bioactivity_type']
                            , value=float(activity_info['value'])
                            , operator=activity_info['operator']
                            , activity_comment=activity_info['activity_comment']
                            , target_confidence=activity_info['target_confidence']
                        )
                        bioassay_data.save()
                    else:
                        self.exceptions.append(Exception('Creating an RDKit molecule failed -- {0}'.format(activity_info['ingredient_cmpd_chemblid'])))
                        self.skipped_compounds.append(compound_data)
        except IntegrityError as ierr:
            self.fatal_exception = ierr
        except Exception as exp:
            if settings.DEBUG:
                raise exp
            else:
                self.fatal_exception = exp