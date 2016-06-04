import math
import pandas
from rdkit.Chem import PandasTools, Descriptors, MolFromSmiles
from sqlalchemy.orm import sessionmaker

import compound_db_utils.database as db
from compound_db_utils import settings


def _gather_columns(table, col_list):
    columns = []
    for col in col_list:
        columns.append(getattr(table.c, col))
    return columns

def fetch_learning_data(
        datasets
        , datasets_cols = ()
        , bioacitivities_cols = (
            'value',
        )
        , compute_descriptors = False
        , create_rdkit_mols = False
        , col_names_map = ()
        , duplicates_handler = None
):
    DB_CONNECTION, TB_COMPOUNDS, TB_DATASETS, TB_BIOACTIVITIES = db.fetch_all()

    session = sessionmaker(bind=DB_CONNECTION)()

    cols = _gather_columns(TB_BIOACTIVITIES, bioacitivities_cols)
    cols.extend(_gather_columns(TB_DATASETS, datasets_cols))
    cols.append(TB_COMPOUNDS.c.smiles)
    query = session.query(
        *cols
    ).join(TB_COMPOUNDS).join(TB_DATASETS)\
        .filter(
        TB_DATASETS.c.unique_id.in_(datasets)
    )

    # make the DB query and export the data to pandas DataFrame object
    data = pandas.read_sql_query(query.selectable, DB_CONNECTION)
    smiles_col_name =  settings.COMPOUNDS_TABLE + '_smiles'
    ic50_col_name = settings.BIOACTIVITIES_TABLE + '_value'

    # remove duplicate values
    if duplicates_handler:
        duplicates = set(data[smiles_col_name][data[smiles_col_name].duplicated()])
        for smiles in duplicates:
            duplicate_ic50s = data[data[smiles_col_name] == smiles][ic50_col_name]
            ret = duplicates_handler(smiles, duplicate_ic50s)
            data = data[data[smiles_col_name] != smiles]
            if type(ret) != bool and ret != False:
                data.update(
                    pandas.DataFrame(
                        [[smiles, ret]]
                        , columns = [smiles_col_name, ic50_col_name]
                    )
                )

    if compute_descriptors:
        desc_list = Descriptors.descList
        try:
            desc_list = [x for x in desc_list if x[0] in compute_descriptors]
        except TypeError:
            for desc_name, function in desc_list:
                values = []
                for smiles in data[smiles_col_name]:
                    mol = MolFromSmiles(smiles)
                    values.append(function(mol))
                data[desc_name] = values


    if create_rdkit_mols:
        PandasTools.AddMoleculeColumnToFrame(
            data
            , smiles_col_name
            , 'rdmol'
        )

    if col_names_map:
        data.rename(columns=col_names_map, inplace=True)

    return data
