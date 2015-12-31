#!/bin/bash

sed -i "s/RAD_GRID\=[0-9]*/RAD_GRID=$1/g" dimensions.f90
sed -i "s/LNG_GRID\=[0-9]*/LNG_GRID=$2/g" dimensions.f90
