

bsub -J split < GoTChA_pipeline_submit_split.sh

bsub -w 'done(split)' -J filt < GoTChA_pipeline_submit_filt.sh

bsub -w 'done(filt)' -J call_mut < GoTChA_pipeline_submit.sh

