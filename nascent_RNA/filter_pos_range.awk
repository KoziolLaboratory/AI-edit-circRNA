BEGIN { OFS="\t" }
NR==FNR {
    ++cnt[$1]
    beg[$1,cnt[$1]] = $2
    end[$1,cnt[$1]] = $3
    next
}
$1 in cnt {
    for (i=1; i<=cnt[$1]; i++) {
        if ( (beg[$1,i] <= $2) && ($2 <= end[$1,i]) ) {
            print $0
        }
    }
}