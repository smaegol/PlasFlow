
from os import environ
from os.path import dirname, join


SCRIPT_PATH = dirname(realpath(sys.argv[0]))
MODELS_PATH = join(SCRIPT_PATH, 'models')
MODELS_PATH = environ.get('PLASFLOW_MODEL_PATH', MODELS_PATH)
