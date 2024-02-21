#!/bin/bash

module purge > /dev/null 2>&1
module load easybuild
module load GCCcore/12.2.0
module load FFmpeg

ffmpeg -i tension_test_30.mp4 -i tension_test_90.mp4 -filter_complex hstack tension_compare.mp4

