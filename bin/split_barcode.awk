#!/usr/bin/awk -f

# make sure the field separator is a tab
# for both input and output
BEGIN {
  FS="\t"
  OFS="\t"
}

# save commented-out rows to print later
/^#/ {
  # make any fields after the 5th blank
  for (i=6; i<=NF; i++) {
    $i = ""
  }
  
  # set output fields to 5
  NF=5
  last_comment = $0
}

# only do stuff for non-commented rows
!/^#/ {
  
  # make seqs uppercase
  $3 = toupper($3)
  $4 = toupper($4)
  $5 = toupper($5)
  
  # if the barcodes look like 'AGCT:',
  # tack the barcode onto the beginning of the forward primer
  # and set the barcodes field to just ':'
  if ($3 ~ /^[AGCT]+:$/) {
    split($3,tags,":")
    $4 = tags[1] $4
    $3 = ":"
  }

  # conversely, if it looks like ':GCTA',
  # tack it on the start of the reverse primer
  # (I think that's how this works)
  if ($3 ~ /^:[AGCT]+$/) {
    split($3,tags,":")
    $5 = tags[2] $5
    $3 = ":"
  }

  # replace I's with N's in the primer sequences (obitools can't handle I's)
  gsub(/I/,"N",$4)
  gsub(/I/,"N",$5)


  # make any fields after the 5th blank
  for (i=6; i<=NF; i++) {
    $i = ""
  }
  # set output fields to 5
  NF=5

  # print last comment to split file if it hasn't been
  if (!comments[$1]) {
    print last_comment > "bc/" $1 "---" parent "_barcode.tsv"
    comments[$1]++
  }
  
  # print into a file named after the first column and some other optional value
  print > "bc/" $1 "---" parent "_barcode.tsv";
}