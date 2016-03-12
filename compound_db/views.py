from django.contrib import messages
from django.core.urlresolvers import resolve, reverse
from django.shortcuts import render, redirect

from compound_db.forms import AddMolForm
from rdkit import Chem
import compound_db.models as models

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
