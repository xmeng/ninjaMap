NAME=ninjamap
TAG:=$(shell /bin/date +%Y%m%d%H%M%S)
REGISTRY=sunitjain/$(NAME)

all: build push run

build:
	docker image build -t $(NAME):$(TAG) -t $(NAME):latest -f Dockerfile .
	docker image tag $(NAME):$(TAG) $(REGISTRY):$(TAG)
	
push:	
	docker image push $(REGISTRY):$(TAG)
	echo "sunitjain/$(NAME):$(TAG)" > LATEST
	
run:
	docker container run --rm $(NAME):$(TAG)
