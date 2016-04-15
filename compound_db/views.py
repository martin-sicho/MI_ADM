import json
import pickle

from django.contrib import messages
from django.core import serializers
from django.core.urlresolvers import resolve, reverse
from django.http import HttpResponse, Http404
from django.shortcuts import render, redirect
from django.template.loader import render_to_string
from rdkit import Chem
from rdkit.Chem import Draw

import compound_db.models as models
from compound_db.ChEMBLImporter import ChEMBLImporter
from compound_db.forms import AddMolForm, ImportChEMBLMols
from compound_db.utils import parse_filters, search_for_term

JSON_MIME_TYPE = 'application/json'

def autocomplete_api(req):
    if req.is_ajax():
        q = req.GET.get('term', '')
        compounds = search_for_term(q)
        results = []
        for idx, compound in enumerate(compounds):
            drug_json = dict()
            drug_json['id'] = compound.unique_id
            drug_json['label'] = compound.smiles
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
            compounds = search_for_term(req.POST.get('basic_query'))
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
            importer = ChEMBLImporter(
                form.cleaned_data['target']
                , form.cleaned_data['description']
                , {
                    'units' : parse_filters(form, 'units')
                    , 'bioactivity_type' : parse_filters(form, 'bioactivity_types')
                    , 'operator' : parse_filters(form, 'operators')
                    , 'value' : (form.cleaned_data['activity_value_threshold'], form.cleaned_data['activity_value_operator'])
                }
            )

            importer.save_data()

            if not importer.fatal_exception:
                return redirect(reverse('compound_db:home'))
            else:
                pickle.dump(importer, open('failed_import.pickle', 'wb'))
                messages.error(req, 'A fatal exception has occured during the save. Import was cancelled.')

    else:
        form = ImportChEMBLMols()
    return render(req, template_name="compound_db/add_chembl_mols.html", context={
        'form' : form
        , 'page_settings' : {
            'navbar_active' : resolve(req.path_info).url_name
            , 'active_title' : 'Import ChEMBL Dataset'
        }
    })
