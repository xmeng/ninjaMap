experiment=141106
sample_sheet="/Users/xdmeng/Downloads/assembly_scripts_handuo/141106.txt"

: <<'END'
for mapgroup in `cat $sample_sheet | sed 1d | awk '{print $3}' | sort | uniq`; do
    echo "MAPGROUP: ${mapgroup} "
	  for coassembly_group in `$sample_sheet | sed 1d | awk '{print $2}' | sort | uniq`; do
		  #group_name=`grep "${coassembly}" /Users/xdmeng/Downloads/assembly_scripts_handuo/DTH161109_grouping.txt | awk '{print $3}' | uniq `
	    #echo " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" DTH161109_grouping
			sampleList=()
			for SAMPLE in `egrep -w "${coassembly_group}" $sample_sheet | egrep -w "${mapgroup}" | awk '{print $1}'`; do
	        #echo " $SAMPLE  ${coassembly} "
					sampleList+=($SAMPLE) 	#Append sample into array ${group_name}
	    done

			if  [ ${#sampleList[@]} -gt 0 ]
			then

			   echo "export sampleList=${sampleList[@]} scripts/run_megahit_mapping_binning.sh  $experiment ${coassembly_group}"    echo "GROUP SIZE: ${#sampleList[@]}"
			fi
			#echo "____________________________________"

	done

done
END

for mapgroup in `cat $sample_sheet | sed 1d | awk '{print $3}' | sort | uniq`; do
    #echo "*****************"
    #echo "MAPGROUP: ${mapgroup} "
    #echo "***************** "
	  for coassembly_group in `egrep -w "${mapgroup}" $sample_sheet | sed 1d | awk '{print $2}' | sort | uniq`; do
		  #group_name=`grep "${coassembly}" /Users/xdmeng/Downloads/assembly_scripts_handuo/DTH161109_grouping.txt | awk '{print $3}' | uniq `
	    #echo " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" DTH161109_grouping
  			for SAMPLE in `egrep -w "${mapgroup}" $sample_sheet | awk '{print $1}'`; do
  	        #echo "${mapgroup} $SAMPLE  ${coassembly_group} "
            echo "export sample=$SAMPLE; export experiment=$experiment; export group=${coassembly_group}; sbatch /home/users/xdmeng/scripts/batch_mapping.sh"
  					#export sample=$SAMPLE; export experiment=$experiment; export group=${coassembly_group}; sbatch /home/users/xdmeng/scripts/batch_mapping.sh
            echo "sleep 3"
  	    done
    done
done
