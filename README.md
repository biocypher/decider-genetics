# DECIDER genetics knowledge graph

This pipeline creates a BioCypher knowledge graph from synthetic data after the
example of the [DECIDER](https://deciderproject.eu) cohort. This knowledge graph
is then connected to a BioChatter instance that allows querying the graph and
other information (papers in a vector database, information from web APIs such
as OncoKB) via natural language.

## 🐳 Run using Docker

> [!IMPORTANT]
> You need to have Docker installed on your machine to run the following commands. Please go to [Docker](https://docs.docker.com/get-docker/) for instructions.

```{bash}
git clone https://github.com/biocypher/decider-genetics.git
cd decider-genetics
docker compose up -d
```

This will build the KG in the Docker container and start a Neo4j instance as
well as a BioChatter Light web app instance configured to only show the KG query
interface. The Neo4j instance can be accessed at `localhost:7474` and the
BioChatter Light instance at `localhost:8501`.

> [!IMPORTANT]
> For using OpenAI GPT as the language model, you will have to provide your API key as an environment variable (`OPENAI_API_KEY`) in your environment. You can do this using an export command (`export OPENAI_API_KEY=sk-...`) or by adding it to your bash profile; you could also provide it to Docker using an env file. We use GPT-3.5-turbo as the default model.

## Questions

The knowledge graph contains information about patients, genes, variants, drugs,
pathways, and clinical data. The schema of the KG can be seen below the query
interface as a JSON object. You can ask questions in natural language, such as:

- How many patients do we have, and what are their names?

- How many clinically significant (CLNSIG = Pathogenic or Likely_pathogenic)
variants does each patient have?

- Does patient1 have a sequence variant in a gene that is druggable with
evidence level "1"? Which drug? Return unique values.

These are only few of infinitely many possible questions, and some may not
result in a valid query. The BioChatter Light interface allows manual
modification and rerunning of the query for prototyping and debugging.

## ⚙️ Local Installation
```{bash}
git clone https://github.com/biocypher/decider-genetics.git
cd decider-genetics
poetry install
poetry run python create_knowledge_graph.py
```

## 🐳 Docker configuration

This repo also contains a `docker compose` workflow to create the example
database using BioCypher and load it into a dockerised Neo4j instance
automatically. To run it, simply execute `docker compose up -d` in the root 
directory of the project. This will start up a single (detached) docker
container with a Neo4j instance that contains the knowledge graph built by
BioCypher as the DB `neo4j` (the default DB), which you can connect to and
browse at localhost:7474. Authentication is deactivated by default and can be
modified in the `docker_variables.env` file (in which case you need to provide
the .env file to the deploy stage of the `docker-compose.yml`).

Regarding the BioCypher build procedure, the `biocypher_docker_config.yaml` file
is used instead of the `biocypher_config.yaml` (configured in
`scripts/build.sh`). Everything else is the same as in the local setup. The
first container (`build`) installs and runs the BioCypher pipeline, the second
container (`import`) installs Neo4j and runs the import, and the third container
(`deploy`) deploys the Neo4j instance on localhost. The files are shared using a
Docker Volume. This three-stage setup strictly is not necessary for the mounting
of a read-write instance of Neo4j, but is required if the purpose is to provide
a read-only instance (e.g. for a web app) that is updated regularly; for an
example, see the [meta graph
repository](https://github.com/biocypher/meta-graph). The read-only setting is
configured in the `docker-compose.yml` file
(`NEO4J_dbms_databases_default__to__read__only: "false"`) and is deactivated by
default.
