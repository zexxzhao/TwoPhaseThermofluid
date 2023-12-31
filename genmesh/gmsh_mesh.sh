#!/bin/bash



if [ $# -eq 2 ]; then
  geofile=$1
  scale=$2
else
  echo "Wrong input"
  echo "Usage: generate_gmsh_mesh <geo-file> <scale>" 
  exit
fi 
 
rm -f gmsh.mesh log
time gmsh $geofile -v 10 -3 -o gmsh.mesh  -format mesh  -clscale $scale | tee log 

echo "scale=$scale" >> log

if [ -f gmsh.mesh ]; then
  nodes=`head -5 gmsh.mesh | tail -1`
  dir=mesh`echo $nodes`
  mkdir -p $dir
  cd $dir
  mv ../gmsh.mesh .
  mv ../log gmsh.log

  grep ille gmsh.log
  echo =============================
  grep ille gmsh.log | grep -v no
  echo =============================
  echo $nodes nodes 
  echo =============================
fi
