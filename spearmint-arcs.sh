#!/bin/sh
set -eux

mkdir -p spearmint-arcs.mongodb && nohup mongod --logpath spearmint-arcs.mongodb/mongod.log --dbpath spearmint-arcs.mongodb >spearmint-arcs.mongodb/nohup.out &
nohup spearmint --config=spearmint-arcs.json . >spearmint-arcs.log &
