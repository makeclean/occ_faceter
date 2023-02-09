# useful commands

(build commands are to be run from this directory)

```shell
# main command to make stripped down image
docker build --pull --tag=occ_facet_geom .

# when debugging it can help to get the "inner" image
docker build --tag=inner --target inner .

# use the image with data in the current directory appearing in /data
docker run -v "$PWD:/data" --rm -it occ_facet_geom

# if you're debugging and want the current filesystem saved for later:
#  ps shows currently running containers
docker ps
#  commit is used to save into new container, where container_name comes from the above line
#  and output will be written to occf_debug
docker commit container_name occf_debug
# get back into it later
docker run --rm it occf_debug
```
