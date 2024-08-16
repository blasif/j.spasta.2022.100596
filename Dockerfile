FROM rocker/geospatial:4.1.0

RUN apt-get update && \
	apt purge texlive* -y && \
	rm -rf /usr/local/texlive/* && \
	rm -rf ~/.texlive* && \
	rm -rf /usr/local/share/texmf && \
	rm -rf /var/lib/texmf && \
	rm -rf /etc/texmf && \
	apt-get remove tex-common --purge -y && \
	rm -rf ~/.texlive && \
	find -L /usr/local/bin/ -lname /usr/local/texlive/*/bin/* | xargs -r rm

ENV DISABLE_AUTH=true
