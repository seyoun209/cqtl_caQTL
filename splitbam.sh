for i in {2,4,5}
do
	echo $i
	sbatch qc_atac.sbatch $i
done
