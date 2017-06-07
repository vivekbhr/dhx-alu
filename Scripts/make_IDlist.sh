
## Human
## Make ID list for cluster annotation (using repeatmasker output), in repeatX pipeline
type=${1}
folder=${2}

if [ type == "human" ]
then
	for seq in AluS AluJ AluY
	do
	nlines=$(grep "${seq}" ${folder}/RM-reads.fas.out | wc -l)
	echo -e "\n"'>'${seq}"\t"${seq}"\t"${nlines} >> ${folder}/IDlist_annotated.txt
	grep "${seq}" ${folder}/RM-reads.fas.out | awk '{print $5}' | tr [:space:] "\t" >> ${folder}/IDlist_annotated.txt
	done
elif [ type == "mouse" ]
then
	for seq in B1_Mur B1_Mus B1_Rn B1_Mm
	do
	nlines=$(grep "${seq}" ${folder}/RM-reads.fas.out | wc -l)
	echo -e "\n"'>'${seq}"\t"${seq}"\t"${nlines} >> ${folder}/IDlist_annotated.txt
	grep "${seq}" ${folder}/RM-reads.fas.out | awk '{print $5}' | tr [:space:] "\t" >> ${folder}/IDlist_annotated.txt
	done
elif [ type == "fly" ]
then
	for seq in Copia
	do
	nlines=$(grep "${seq}" ${folder}/RM-reads.fas.out | wc -l)
	echo -e "\n"'>'${seq}"\t"${seq}"\t"${nlines} >> ${folder}/IDlist_annotated.txt
	grep "${seq}" ${folder}/RM-reads.fas.out | awk '{print $5}' | tr [:space:] "\t" >> ${folder}/IDlist_annotated.txt
done
fi
