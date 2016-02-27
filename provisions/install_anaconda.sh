function install_anaconda {
  INSTALLER=$1
  INSTALLATION_DIR=$2

  NEW_PATH="$INSTALLATION_DIR/bin:$PATH"
  if [ -s $INSTALLER ]; then
    chmod +x $INSTALLER
    if $INSTALLER -f -b -p $INSTALLATION_DIR; then
      echo "Anaconda installed in $INSTALLATION_DIR"
    else
      printf '%s\n' 'Anaconda installation failed!' >&2
      return 1
    fi
    echo "Appending Anaconda to path..."
    touch $HOME/.bashrc
    echo -e "\n#For anaconda\nPATH=$NEW_PATH" >> $HOME/.bashrc
  else
    echo "ERROR: Anaconda installer not found"
    return 1
  fi
}

install_anaconda $1 $2