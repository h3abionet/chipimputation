#!/usr/bin/make -f

# generate ca and make keys from itx

SWARM_MACHINES:=manager0 manager1 consul0 \
	$(patsubst %,node%,$(shell seq -w 0 99))

all: $(patsubst %,%.cert,$(SWARM_MACHINES)) $(patsubst %,%-priv-key.pem,$(SWARM_MACHINES))

%-priv-key.pem:
	openssl genrsa -out $@ 2048

## this is the maximum number of days that keys can be valid for
## before we hit the unix epoch
MAX_DAYS:=$(shell echo $$(( ( `date -d "Tuesday, 18 January 2038" +%s` - `date +%s` ) / 60 / 60 /24 )))

## generate the ca.pem and a database/config for it
ca.pem: ca-priv-key.pem
	if [ ! -e database ]; then \
		cp /usr/lib/ssl/openssl.cnf config; \
		touch database database.attr; \
		perl -pi -e 's/(database|serial)\s*=.+/$$1=$$1/' config; \
		perl -e 'print unpack(q(H*),pack(q(NN),time,$$$$)),qq(\n)' > serial; \
	fi;
	echo -e ".\n.\n.\n.\n.\nImputationSwarm\n.\n" | \
		openssl req -config ./config -new -key $< -x509 -days $(MAX_DAYS) -out $@;

# generate the Certificate Signing Request
%.req: %-priv-key.pem
	openssl req -subj "/CN=$*" -new -key $< -days $(MAX_DAYS) -out $@

%.cert: %.req ca-priv-key.pem ca.pem
	openssl ca -config config -policy policy_anything -keyfile $(wordlist 2,2,$^) \
		-cert $(wordlist 3,3,$^) -out $@ -outdir . -notext -days $(MAX_DAYS) \
		-batch -infiles $<
	chmod a+r $@


clean:
	rm -f *.cert *.pem
