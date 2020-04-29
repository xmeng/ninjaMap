NAME=ninjamap
TAG:=$(shell /bin/date +%Y%m%d%H%M%S)
REGISTRY=sunitjain/$(NAME)

all: dev
master: dev m_build m_push m_run
dev: d_build d_push d_run

d_build:
	docker image build --pull --no-cache -t $(NAME):$(TAG) -f Dockerfile .
	docker image tag $(NAME):$(TAG) $(REGISTRY):$(TAG)
	docker container run --rm $(NAME):$(TAG) samtools --help

d_push:	
	docker image push $(REGISTRY):$(TAG)
	echo "sunitjain/$(NAME):$(TAG)" > LATEST
	
d_run:
	docker container run --rm $(NAME):$(TAG)

m_build:
	# docker image build -t $(NAME):latest -f Dockerfile .
	docker image tag $(NAME):$(TAG) $(REGISTRY):latest
	
m_push:	
	docker image push $(REGISTRY):latest
	
m_run:
	docker container run --rm $(REGISTRY):latest