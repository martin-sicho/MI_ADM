import json

from django.contrib import messages
from django.core.urlresolvers import resolve, reverse
from django.http import HttpResponse
from django.shortcuts import render, redirect

from compound_db.forms import AddMolForm
from rdkit import Chem
import compound_db.models as models

def search_api(req):
    if req.is_ajax():
        q = req.GET.get('term', '')
        compounds = models.Compound.objects.filter(unique_id__icontains = q)[:20]
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
    mimetype = 'application/json'
    return HttpResponse(data, mimetype)

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
        form = AddMolForm()
    return render(req, template_name="compound_db/add_mol.html", context={
        'form' : form
        , 'page_settings' : {
            'navbar_active' : resolve(req.path_info).url_name
            , 'active_title' : 'Add Molecule'
        }
    })
