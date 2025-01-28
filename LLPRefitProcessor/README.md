Based on RefitProcessor in MarlinTrkProcessors

### Build procedure

cmake -C $ILCSOFT/ILCSoft.cmake ..
make install

### Include in used processors ###

export MARLIN_DLL=./lib/libLLPRefitProcessor.so
