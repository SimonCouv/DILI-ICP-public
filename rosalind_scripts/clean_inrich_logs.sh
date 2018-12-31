cd "/mnt/lustre/users/k1893262/DILI-ICP/DILI-ICP/rosalind_logs/"
logfile="inrich.o4382747"

grep -vP "Reading \d+ map positions" "$logfile" | grep -vP "\d+ SNP counts assigned" | grep -vP "Precomputing acceptable positions \d+" | grep -vP "\d+ first-pass permutations" | grep -vP "\d+ second-pass permutations" > inrich_cleaned.log