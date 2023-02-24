singularity pull vina_1.2.3.0.sif docker://tiffhuff/autodock_vina:1.2.3.0

flex=${flexible_sidechains}
if
     [ -z $flex ]
then
     flex='Empty'
fi
    
MV2_ENABLE_AFFINITY=0 ibrun singularity exec vina_1.2.3.0.sif python3 autodock.py \
     -r ${receptor} \
     -c "${center_x},${center_y},${center_z}" \
     -s "${size_x},${size_y},${size_z}" \
     -m ${forcefield} \
     -d ${docking} \
     -ll ${library} \
     -n ${top_n_scores} \
     -f ${flex}