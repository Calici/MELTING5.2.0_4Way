FROM debian:bullseye-slim AS git-cloner

RUN apt-get -y update
RUN apt-get -y install git


RUN git clone https://github.com/Calici/MELTING5.2.0_4Way


FROM openjdk:21-ea-17-slim-bullseye AS java-exec


RUN apt-get -y update
RUN apt-get -y install python3

COPY --from=git-cloner /MELTING5.2.0_4Way /MELTING5.2.0_4Way
WORKDIR /MELTING5.2.0_4Way/executable

RUN chmod +x melting


ENTRYPOINT [ "python3", "4wayJunction.py" ]
