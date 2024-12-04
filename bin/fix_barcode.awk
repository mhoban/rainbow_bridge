#!/usr/bin/awk -f
BEGIN { FS="\t"; OFS="\t" }
/^#/ {print} 
!/^#/ {gsub(/[Ii]/,"N",$4);gsub(/[Ii]/,"N",$5);print}