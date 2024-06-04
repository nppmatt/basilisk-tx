#!/bin/bash

BIN_DIR=bin
options=( "$OUT_DIR"/* )
select option in "${options[@]}" "cancel"; do
	case $option in
		"cancel")
			break
			;;
		*)
			echo "Scheduling $option"
			sbatch run_job.sh "$option"
			break
			;;
	esac
done

echo "Scheduling done."

