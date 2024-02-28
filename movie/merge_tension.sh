#!/bin/bash

module purge > /dev/null 2>&1
module load easybuild
module load GCCcore/12.2.0
module load FFmpeg

ffmpeg -i tension_test_15.mp4 -i tension_test_100.mp4 -filter_complex hstack compare_15_100.mp4

