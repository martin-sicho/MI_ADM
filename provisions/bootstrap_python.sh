function install_python_packages {
  CONDA_BIN=$1/bin
  PATH="$CONDA_BIN:$PATH"
  echo $PATH
  conda install -y -f -c https://conda.anaconda.org/rdkit django rdkit bootstrap3 && pip install django-bootstrap3
}

install_python_packages $1