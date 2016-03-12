from django import forms

class AddMolForm(forms.Form):

    mol_file = forms.CharField(
        widget=forms.HiddenInput()
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