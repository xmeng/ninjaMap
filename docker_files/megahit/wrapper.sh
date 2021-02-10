myexperiment="pilot"
sample_sheet="/Users/xdmeng/Downloads/assembly_scripts_handuo/pilot-3col.txt"

for mapgroup in `cat $sample_sheet | sed 1d | awk '{print $3}' | sort | uniq`; do
    echo "*****************"
    echo "MAPGROUP: ${mapgroup} "
    echo "***************** "
	  for coassembly_group in `cat $sample_sheet | sed 1d | awk '{print $2}' | sort | uniq`; do
		  #group_name=`grep "${coassembly}" /Users/xdmeng/Downloads/assembly_scripts_handuo/DTH161109_grouping.txt | awk '{print $3}' | uniq `
	    #echo " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" DTH161109_grouping
			sampleList=()
			for SAMPLE in `egrep -w "${coassembly_group}" $sample_sheet | egrep -w "${mapgroup}" | awk '{print $1}'`; do
	        #echo " $SAMPLE  ${coassembly} "
					sampleList+=($SAMPLE) 	#Append sample into array ${group_name}
	    done

			if  [ ${#sampleList[@]} -gt 0 ]
			then
         echo "GROUP SIZE: ${#sampleList[@]}"
         mysamples=$( IFS=:; printf '%s' "${sampleList[*]}" )
         #echo "export samples=$mysamples; export experiment=$myexperiment; export group=${coassembly_group}; sbatch /home/users/xdmeng/scripts/batch_assembly.sh"
			   #echo "export sampleList=(${sampleList[@]}) scripts/run_megahit_mapping_binning.sh $experiment ${coassembly_group}"
         #echo "---------------------------------"
			fi
			#echo "____________________________________"

	done

done
