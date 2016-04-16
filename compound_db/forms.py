import os

from chembl_webresource_client import TargetResource
from django import forms
from django.conf import settings

from compound_db.utils import load_filter_choices


class AddMolForm(forms.Form):

    mol_file = forms.CharField(
        widget=forms.HiddenInput(attrs={'class': "mol_file_field"})
    )
    description = forms.CharField(
        widget = forms.Textarea(
            attrs={
                'placeholder' : 'A short text describing this molecule.'
                , 'class' : 'materialize-textarea'
            }
        )
        , required=False
    )

class ImportChEMBLMols(forms.Form):

    target = forms.CharField(
        required=True
        , widget=forms.TextInput(
            attrs={
                'placeholder' : 'ChEMBL target ID (e.g. CHEMBL206)'
            }
        )
    )

    description = forms.CharField(
        widget = forms.Textarea(
            attrs={
                'placeholder' : 'A short text describing the purpose of this data.'
                , 'class' : 'materialize-textarea'
            }
        )
        , required=False
    )

    units = forms.MultipleChoiceField(
        widget=forms.SelectMultiple
        , required=False
        , choices=load_filter_choices(os.path.join(settings.BASE_DIR, 'compound_db/units.txt'))
    )
    units_custom = forms.CharField(
        required=False
        , widget=forms.TextInput(
            attrs={
                'placeholder' : 'Other units to include in the final set (comma seperated).'
            }
        )
        , label='Custom Units'
    )

    bioactivity_types = forms.MultipleChoiceField(
        widget=forms.SelectMultiple
        , required=False
        , choices=load_filter_choices(os.path.join(settings.BASE_DIR, 'compound_db/bioactivity_types.txt'))
    )
    bioactivity_types_custom = forms.CharField(
        required=False
        , widget=forms.TextInput(
            attrs={
                'placeholder' : 'Other bioactivity types to include in the final set (comma seperated).'
            }
        )
        , label='Custom Bioactivity Types'
    )

    operators = forms.MultipleChoiceField(
        widget=forms.SelectMultiple
        , required=False
        , choices=load_filter_choices(os.path.join(settings.BASE_DIR, 'compound_db/operator_types.txt'))
    )
    operators_custom = forms.CharField(
        required=False
        , widget=forms.TextInput(
            attrs={
                'placeholder' : 'Other operators to include in the final set (comma seperated).'
            }
        )
        , label='Custom Operators'
    )

    activity_value_threshold = forms.FloatField(
        required=False
        , label="Value"
        , widget=forms.TextInput(
            attrs={
                'placeholder' : 'Any number'
            }
        )
    )
    activity_value_operator = forms.ChoiceField(
        choices=load_filter_choices(
            os.path.join(settings.BASE_DIR, 'compound_db/threshold_operators.txt')
            , prepend=None
        )
        , widget=forms.RadioSelect()
        , label="Operator"
        , initial="<="
    )

    activity_comments = forms.CharField(
        required=False
        , widget=forms.TextInput(
            attrs={
                'placeholder' : 'Desired activity comments (comma seperated).'
            }
        )
        , label='Activity Comment'
    )

    def clean(self):
        cleaned_data = super(ImportChEMBLMols, self).clean()
        if 'target' in cleaned_data:
            chembl_target_id = cleaned_data['target']

            targets = TargetResource()
            if not targets.status():
                raise forms.ValidationError(
                        "The ChEMBL web service is currently unavailable. "
                        "Try submitting your request later."
                    )

            target = targets.get(chembl_target_id)
            if type(target) != dict or not target:
                self.add_error('target',
                               "This target record does not seem to exist on ChEMBL. "
                               "Make sure you are submitting a correct ID.")

        return cleaned_data