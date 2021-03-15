#!/bin/bash -x
## -x : expands variables and prints commands as they are called

## TODO:
## - specification: after -- (empty flag) use params not used for this script -> i.e. reserved for tool not part of this image
## - make flags for each tool -> if set, create index for said tool
## - (Salmon) replace Salmon call
## - (R)remove organism and taxonomy id from R index?
## - remove hisat2_index.sh
## - (kallisto) replace kallisto absolute path with alias
## - (kallisto/Salmon) resolve kallisto/Salmon dependency
## - (kallisto) update kallisto
## - (STAR) replace overhang with readLength param provided from automated config
## - echo mapped variables before image call
## - print mapping /home/blabla/index mapped to /home/data/index


# saner programming env: these switches turn some bugs into errors
set -o errexit -o pipefail -o noclobber -o nounset

# -allow a command to fail with !’s side effect on errexit
# -use return value from ${PIPESTATUS[0]}, because ! hosed $?
! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo 'I’m sorry, `getopt --test` failed in this environment.'
    exit 1
fi

OPTIONS=
LONGOPTS=gtf:,fasta:,organism:,taxid:,nthread:,hisat2,star,kallisto,salmon,r,dexseq,index:,empires,log:

# -regarding ! and PIPESTATUS see above
# -temporarily store output to be able to check for errors
# -activate quoting/enhanced mode (e.g. by writing out “--options”)
# -pass arguments only via   -- "$@"   to separate them correctly
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi
# read getopt’s output this way to handle the quoting right:
eval set -- "$PARSED"

gtf=- fasta=- organism=- taxid=- nthread=4
hisat2=n star=n kallisto=n salmon=n r=n dexseq=n
# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
        --gtf)
            gtf="$2"
            shift 2
            ;;
        --fasta)
            fasta="$2"
            shift 2
            ;;
        --organism)
            organism="$2"
            shift 2
            ;;
        --taxid)
            taxid="$2"
            shift 2
            ;;
        --nthread)
            nthread="$2"
            shift 2
            ;;
        --hisat2)
            hisat2=y
            shift
            ;;
        --star)
            star=y
            shift
            ;;
        --kallisto)
            kallisto=y
            shift
            ;;
        --salmon)
            salmon=y
            shift
            ;;
        --dexseq)
            dexseq=y
            shift
            ;;
		--r)
            r=y
            shift
            ;;
        --)
            shift
            break
            ;;
		--log)
			log="$2"
			shift 2
			;;
        *)
            shift
            ;;
    esac
done

# handle non-option arguments
if [[ $# -ne 0 ]]; then
    echo "$0: empty flag detected, is this intentional?!"
    #exit 4
fi


outdir='/home/data/indices'
echo 'organism:'$'\t'$organism
echo 'taxid:'$'\t'$taxid
echo 'gtf:'$'\t'$gtf
echo 'fasta:'$'\t'$fasta

watch pidstat -du -hl '>>' $log/star-$(date +%s).pidstat & wid=$!

## STAR
if [[ "$star" = "y" ]]; then
    echo $'\n'"[INFO] [generate_indices.sh] [STAR] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Start: generate STAR index..."
    overhang=99
    mkdir -p $outdir/STAR/$overhang
    /usr/bin/time -v /home/software/STAR/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --runThreadN $nthread \
    --genomeDir $outdir/STAR/ --genomeFastaFiles $fasta \
    --sjdbGTFfile $gtf --sjdbOverhang $overhang
    echo $'\n'"[INFO] [generate_indices.sh] [STAR] ["`date "+%Y/%m/%d-%H:%M:%S"`"] End: generate STAR index..."
fi

kill -15 $wid