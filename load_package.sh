#!/bin/bash

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "This script should be sourced, not executed. Use 'source ${0}'"
    exit 1
fi

echo "Select an option to activate:"
echo "1: Activate ROOT"
echo "2: Activate SVFIT"
echo "3: Activate OPTION3"
echo "4: Activate OPTION4"

read -p "Enter the option number: " option

case $option in
    1)
        echo "Load ROOT..."
        source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/thisroot.sh
        ;;
    2)
        echo "Load SVFIT..."
        LD_LIBRARY_PATH=/usr/local/lib
        CPPYY_BACKEND_LIBRARY=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/libcppyy_backend3_9.so
        source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/bin/thisroot.sh
        export LIBRARY_PATH=$LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/TauAnalysis/ClassicSVfit/lib
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.30.04/x86_64-centosstream9-gcc113-opt/lib/
        python -c "import TauAnalysis.ClassicSVfit.wrapper.pybind_wrapper" && echo "Module is working" || echo "Module import failed"
        ;;
    3)
        echo "Activating OPTION3..."
        ;;
    4)
        echo "Activating OPTION4..."
        ;;
    *)
        echo "Invalid option selected."
        ;;
esac
