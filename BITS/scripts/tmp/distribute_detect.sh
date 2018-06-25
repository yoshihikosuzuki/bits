#!/bin/bash

SPLIT_NUM=$1
DALIGNER_NCORE=$2
QSUB_QUEUE=$3

pyfasta split -n ${SPLIT_NUM} preads4falcon.fasta

for i in $(seq -f %03g 0 $(expr ${SPLIT_NUM} - 1)); do
	echo -e "#!/bin/bash\nsource ~/work2/pbmga/pbmga_env/bin/activate\npython ~/work2/pbmga/src/detect_multiplexed_sequences.py --out_file preads.adaptor_unremoved.${i} preads4falcon.${i}.fasta" > qsub_detect.${i}.sh
	qsub -N detect${i} -o sge.log -j y -pe smp ${DALIGNER_NCORE} -q ${QSUB_QUEUE} -S /bin/bash -cwd qsub_detect.${i}.sh
done

echo -e "#!/bin/bash\ncat preads.adaptor_unremoved.* > preads.adaptor_unremoved" > qsub_gather.sh
JID_LIST="detect000"
for i in $(seq -f %03g 1 $(expr ${SPLIT_NUM} - 1)); do
	JID_LIST="${JID_LIST},detect${i}"
done
qsub -N gather -o sge.log -j y -q ${QSUB_QUEUE} -S /bin/bash -cwd -sync y -hold_jid ${JID_LIST} qsub_gather.sh

rm preads.adaptor_unremoved.*
rm preads4falcon.*.fasta
rm qsub_*
