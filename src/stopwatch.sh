#!/bin/bash

START=$SECONDS

$1

END=$SECONDS

echo -e "\n"$((END - START))
