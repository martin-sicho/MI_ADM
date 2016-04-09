import json

from chembl_webresource_client import TargetResource, CompoundResource
from django.contrib import messages
from django.core import serializers
from django.core.urlresolvers import resolve, reverse
from django.db import transaction, IntegrityError
from django.http import HttpResponse, Http404
from django.shortcuts import render, redirect
from django.template.loader import render_to_string
from rdkit.Chem import Draw

from compound_db.forms import AddMolForm, ImportChEMBLMols
from rdkit import Chem
import compound_db.models as models
from compound_db.utils import is_number

JSON_MIME_TYPE = 'application/json'

CHEMBL_TARGETS_RESOURCE = TargetResource()
CHEMBL_COMPOUNDS_RESOURCE = CompoundResource()

def autocomplete_api(req):
    if req.is_ajax():
        q = req.GET.get('term', '')
        compounds = models.Compound.objects.filter(unique_id__icontains = q)[:5]
        results = []
        for idx, compound in enumerate(compounds):
            drug_json = dict()
            drug_json['id'] = compound.unique_id
            drug_json['label'] = compound.unique_id
            drug_json['value'] = compound.unique_id
            results.append(drug_json)
        data = json.dumps(results)
    else:
        data = 'fail'
    return HttpResponse(data, JSON_MIME_TYPE)

def search_api(req):
    if req.method == 'POST' and req.is_ajax():
        data = json.dumps(req.POST)
        if 'basic_query' in req.POST:
            compounds = models.Compound.objects.filter(unique_id__icontains = req.POST.get('basic_query'))[:10]
            data = json.loads(serializers.serialize("json", compounds))
            for item in data:
                del item['model']
            ret = dict()
            ret['data'] = data
            ret['html'] = render_to_string("compound_db/compound_table.html", request=req, context={
                'compounds' : compounds
                , 'column_names' : ['ID', 'SMILES', 'Weight', 'Image']
            })
            return HttpResponse(json.dumps(ret), JSON_MIME_TYPE)
        else:
            return Http404('wrong query')
    else:
        raise Http404

def molimages_api(req, unique_id):
    compound = models.Compound.objects.get(unique_id=unique_id)
    mol = Chem.MolFromSmiles(compound.smiles)
    img = Draw.MolToImage(mol, size=(50,50))
    response = HttpResponse(content_type="image/svg")
    img.save(response, "PNG")
    return response

def home(req):
    return render(req, template_name="compound_db/home.html", context={
        'page_settings' : {
            'navbar_active' : resolve(req.path_info).url_name
            , 'active_title' : 'Search'
        }
    })

def add_mol(req):
    form = None
    if req.method == 'POST':
        form = AddMolForm(req.POST)
        if form.is_valid():
            description = form.cleaned_data['description']
            mol = Chem.MolFromMolBlock(form.cleaned_data['mol_file'])
            if mol:
                new_mol = models.Compound(mol, description)
                try:
                    new_mol.save()
                    messages.success(req, 'Molecule saved successfully. You should now find it in the database.')
                    return redirect(reverse('compound_db:home'))
                except models.Compound.MoleculeAlreadyExists as exp:
                    messages.error(req, str(exp))
                # except Exception:
                #     messages.error(req, 'Unknown error has occured.')
            else:
                raise Exception('RDKit failed to read molecule.')
    else:
        form = AddMolForm()
    return render(req, template_name="compound_db/add_mol.html", context={
        'form' : form
        , 'page_settings' : {
            'navbar_active' : resolve(req.path_info).url_name
            , 'active_title' : 'Add Molecule'
        }
    })

def add_chembl_mols(req):
    form = None
    if req.method == 'POST':
        form = ImportChEMBLMols(req.POST)
        if form.is_valid():
            # TODO: test if connection OK

            skipped_compounds = [] # TODO: log this in the database for future review
            with transaction.atomic():
                target = form.cleaned_data['target']
                target_data = models.ChEMBLTargetData(
                    uniprot_accession=target['proteinAccession']
                    , chembl_id=target['chemblId']
                    , organism=target['organism']
                    , preffered_name=target['preferredName']
                    , description=form.cleaned_data['description']
                )
                target_data.save()

                activities = CHEMBL_TARGETS_RESOURCE.bioactivities(target['chemblId'])
                for idx,compound_activity in enumerate(activities):
                    print('Processing {0}/{1}...'.format(idx, len(activities)))
                    compound_data = CHEMBL_COMPOUNDS_RESOURCE.get(compound_activity['ingredient_cmpd_chemblid'])
                    if type(compound_data) != dict \
                            or 'smiles' not in compound_data \
                            or not is_number(compound_activity['value']):
                        skipped_compounds.append(compound_data)
                        continue
                    mol = Chem.MolFromSmiles(compound_data['smiles'])
                    if mol:
                        compound = models.Compound(mol, '')
                        try:
                            with transaction.atomic():
                                compound.save()
                        except models.Compound.MoleculeAlreadyExists:
                            compound = models.Compound.objects.get(inchi_key=compound.inchi_key)
                        except IntegrityError:
                            skipped_compounds.append(compound_data) # FIXME: this results in a loss of activity data if an enantiomer is already in the database
                            continue

                        bioassay_data = models.ChEMBLBioassayData(
                            compound=compound
                            , target_data=target_data
                            , assay_id=compound_activity['assay_chemblid']
                            , ingredient_cmpd_id=compound_activity['ingredient_cmpd_chemblid']
                            , units=compound_activity['units']
                            , bioactivity_type=compound_activity['bioactivity_type']
                            , value=float(compound_activity['value'])
                            , operator=compound_activity['operator']
                            , activity_comment=compound_activity['activity_comment']
                            , target_confidence=compound_activity['target_confidence']
                        )
                        bioassay_data.save()
                    else:
                        skipped_compounds.append(compound_data)

            return redirect(reverse('compound_db:home'))

    else:
        form = ImportChEMBLMols()
    return render(req, template_name="compound_db/add_chembl_mols.html", context={
        'form' : form
        , 'page_settings' : {
            'navbar_active' : resolve(req.path_info).url_name
            , 'active_title' : 'Import ChEMBL Dataset'
        }
    })
