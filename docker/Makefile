DOCKER_IMAGE_NAME := soapdenovo-trans
DOCKER_IMAGE_TAG ?= latest

help:
	@ echo "build :  build $(DOCKER_IMAGE_NAME) binaries with a docker container"
	@ echo "clean :  clean the binaries"

.PHONY: docker-version
docker-version:
	docker version

.PHONY: docker-image-check
docker-image-check:
	@if [ -f docker-image -a -z "$$(docker images -q $(DOCKER_IMAGE_NAME):$(DOCKER_IMAGE_TAG))" ]; then \
		echo "docker-image file exists and docker image does not, removing check file"; \
		rm docker-image; \
	fi

docker-image: Dockerfile
	docker build -t $(DOCKER_IMAGE_NAME):$(DOCKER_IMAGE_TAG) .
	@ touch $@

build: docker-version docker-image-check docker-image
	@echo "Executing 'make.sh' in docker"
	docker run --rm \
		-v $$(pwd)/..:/app \
		-w /app \
		$(DOCKER_IMAGE_NAME):$(DOCKER_IMAGE_TAG) \
		bash make.sh

clean: docker-version docker-image-check docker-image
	@echo "Executing 'clean.sh' in docker"
	docker run --rm \
		-v $$(pwd)/..:/app \
		-w /app \
		$(DOCKER_IMAGE_NAME):$(DOCKER_IMAGE_TAG) \
		bash clean.sh