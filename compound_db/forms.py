from chembl_webresource_client import TargetResource
from django import forms

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
        , choices=(
            ('All', 'All')
            , ('nM', 'nM')
            , ('uM', 'uM')
            , ('percentage', '%')
            , ('Unspecified', 'Unspecified')
        )
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