
DATABASE = {
    'NAME': 'MI_ADM',
    'USER': 'MI_ADM',
    'PASSWORD': 'MI_ADM',
    'HOST': '127.0.0.1',
    'PORT': '5432',
    'ENGINE' : 'psycopg2'
}

TABLE_PREFIX = 'compound_db_'
COMPOUNDS_TABLE = TABLE_PREFIX + 'compound'
DATASETS_TABLE = TABLE_PREFIX + 'chemblbioassaydata'
BIOACTIVITIES_TABLE = TABLE_PREFIX + 'chembltargetdata'