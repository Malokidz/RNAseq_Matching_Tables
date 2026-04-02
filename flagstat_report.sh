echo -e "Sample\tTotal\tMapped\tMapped_%\tProperly_Paired\tProperly_Paired_%" > merged_flagstat.tsv

for f in *_sorted.flagstat.txt
do
    sample=$(basename $f _sorted.flagstat.txt)

    total=$(grep "in total" $f | awk '{print $1}')
    mapped=$(grep "mapped (" $f | head -1 | awk '{print $1}')
    mapped_pct=$(grep "mapped (" $f | head -1 | awk -F '[()%]' '{print $2}')
    proper=$(grep "properly paired" $f | awk '{print $1}')
    proper_pct=$(grep "properly paired" $f | awk -F '[()%]' '{print $2}')

    echo -e "$sample\t$total\t$mapped\t$mapped_pct\t$proper\t$proper_pct" >> merged_flagstat.tsv
done
