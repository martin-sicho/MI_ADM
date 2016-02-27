VAGRANT_HOME=/home/vagrant/
PROVISIONS_DIR=$VAGRANT_HOME/MI_ADM/provisions/
INSTALLATION_DIR=$VAGRANT_HOME/anaconda/
cd $PROVISIONS_DIR

if su -l vagrant -c "$PROVISIONS_DIR/install_anaconda.sh $PROVISIONS_DIR/files/anaconda.sh $INSTALLATION_DIR"; then
  echo "Successfully installed and set up Anaconda."
else 
  printf '%s\n' 'Failed to install Anaconda' >&2
  exit 1
fi

su -l vagrant -c "$PROVISIONS_DIR/bootstrap_python.sh $INSTALLATION_DIR"
