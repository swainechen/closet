#!/bin/tcsh
cd /home/slchen/mail
mv cron cron-`date +%F` && touch cron && bzip2 cron-`date +%F`
