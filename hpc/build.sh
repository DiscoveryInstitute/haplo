#-----------------------------------------------------------------------------
#
#  Compile the code
#
#-----------------------------------------------------------------------------

ROOTPATH=$( cd $(dirname $0)"/.." ; pwd -P)
SCRIPTPATH=$ROOTPATH"/hpc"
COMPILEPATH=$ROOTPATH"/scratch/hpc/default"

mkdir -p $COMPILEPATH
cd $COMPILEPATH

#make -s -f $SCRIPTPATH/Makefile clean
make -s -j $(nproc) -k -f $SCRIPTPATH/Makefile $1 2> err
head -n 25 err
# -s do not echo everything
# -j use multiple procs
# -f name of makefile (default: Makefile) 
# -k says keep going after errors

#-----------------
echo Compiled to: $COMPILEPATH
cd - > /dev/null
#-----------------
