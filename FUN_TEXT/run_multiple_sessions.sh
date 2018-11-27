# loop through all files in /DATA/AZEL/
for file in `ls ../DATA/AZEL/azel_*.txt | xargs -n1 basename`
do

    echo $file   # this just outputs the name of the current azel file
	./radiate $file vlbi.ell -one_epoch_per_obs -cleanup -notrp       # this is the actual command which starts the ray-tracing
	
done





% this writes the duration of the processing to the command window
end=date
duration=$SECONDS
echo "Finished after $(($duration/3600)) hours $(($duration/60 % 60)) minutes and $(($duration % 60)) seconds"

