import traceback

from django.db import DataError

from MI_ADM import settings
from compound_db import models


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def load_filter_choices(file, prepend=('All', 'All')):
    ret = list()
    if prepend:
        ret.append(prepend)
    with open(file, 'r', encoding='utf-8') as infile:
        for line in infile:
            line = line.strip().strip("\"")
            if line:
                ret.append((line, line))
    return tuple(ret)

def parse_filters(form, filter_key):
    filter_options = set(form.cleaned_data[filter_key]) | set(form.cleaned_data[filter_key + '_custom'].split(','))
    filter_options.discard('')
    return filter_options

def search_for_term(term, offset=0, max_items=10):
    queries = [
        lambda : models.Compound.objects.filter(unique_id__icontains = term)
        , lambda : models.Compound.objects.filter(inchi__icontains = term)
        , lambda : models.Compound.objects.filter(inchi_key__icontains = term)
        , lambda : models.Compound.objects.filter(mol__hassubstruct = term)
    ]

    ret = list()
    for query in queries:
        if len(ret) <= max_items:
            try:
                compounds = query()
                if offset < len(compounds):
                    ret.extend(compounds[offset:offset+max_items])
            except DataError as err:
                if settings.DEBUG:
                    try:
                        raise err
                    except DataError:
                        pass
                    traceback.print_exc()
                else:
                    pass
        else:
            break

    return ret[:max_items]