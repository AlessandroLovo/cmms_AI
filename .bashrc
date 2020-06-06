
#------------------------------------------------------------------------
# this is for XCRYSDEN 1.4.1; added by XCRYSDEN installation on
# Tue Mar 24 10:01:25 CET 2009
#------------------------------------------------------------------------
# GIT_DIR="/media/alessandro/OS/Users/ale57/Documents/1.\ universita\'/ANNO\ IV\ \(2019-2020\)/second\ semester/Computational\ Methods\ for\ MS/spiro/ABINITIO/Al100"
# export GIT_DIR

XCRYSDEN_TOPDIR=~/XCrySDen/xcrysden-1.6.2-bin-shared
XCRYSDEN_SCRATCH=~/scratch
QE_TOPDIR=~/Quantum-Espresso/qe-6.5
export XCRYSDEN_TOPDIR XCRYSDEN_SCRATCH QE_TOPDIR
PSEUDO_DIR="$QE_TOPDIR/pseudo/"
export PSEUDO_DIR
PATH="$XCRYSDEN_TOPDIR:$PATH:$XCRYSDEN_TOPDIR/scripts:$XCRYSDEN_TOPDIR/util"
PATH=$PATH:/usr/local/bin
PATH=$PATH:/usr/common/
export PATH
alias pw.x="$QE_TOPDIR/bin/pw.x"
alias pp.x="$QE_TOPDIR/bin/pp.x"
alias plotrho.x="$QE_TOPDIR/bin/plotrho.x"
alias rm='\rm -i'

#------------------------------------------------------------------------
# this is for XCRYSDEN 1.4.1; added by XCRYSDEN installation on
# Mon May 11 13:11:37 CEST 2009
#------------------------------------------------------------------------
# XCRYSDEN_TOPDIR=/usr/common/xc
# XCRYSDEN_SCRATCH=/scratch/franc/xcrys_tmp
# export XCRYSDEN_TOPDIR XCRYSDEN_SCRATCH
PATH="$XCRYSDEN_TOPDIR:$PATH:$XCRYSDEN_TOPDIR/scripts:$XCRYSDEN_TOPDIR/util"

