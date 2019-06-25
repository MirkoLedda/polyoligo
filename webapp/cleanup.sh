#!/bin/bash

find uploads/ -mtime +1 | xargs rm -rf
