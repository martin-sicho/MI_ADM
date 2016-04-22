
DATABASE = { # only postgresql supported so far
    'NAME': 'MI_ADM',
    'USER': 'MI_ADM',
    'PASSWORD': 'MI_ADM',
    'HOST': '127.0.0.1',
    'PORT': '5432',
    'ENGINE' : 'psycopg2'
}

TABLE_PREFIX = 'compound_db_'

COMPOUNDS_TABLE = TABLE_PREFIX + 'compound'
COMPOUNDS_TABLE_COLS = [
    'id'
    , 'unique_id'
    , 'smiles'
]

DATASETS_TABLE = TABLE_PREFIX + 'chembltargetdata'
DATASETS_TABLE_COLS = [
    'id'
    , 'unique_id'
    , 'chembl_id'
    , 'organism'
    , 'preffered_name'
    , 'description'
]

BIOACTIVITIES_TABLE = TABLE_PREFIX + 'chemblbioassaydata'
BIOACTIVITIES_TABLE_COLS = [
    'id'
    , 'compound_id'
    , 'target_data_id'
    , 'assay_id'
    , 'ingredient_cmpd_id'
    , 'units'
    , 'bioactivity_type'
    , 'value'
    , 'operator'
    , 'activity_comment'
    , 'target_confidence'
]