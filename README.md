# DECIDER genetics knowledge graph and decision support

This pipeline creates a BioCypher knowledge graph from synthetic data after the
example of the [DECIDER](https://deciderproject.eu) cohort. This knowledge graph
is then connected to a BioChatter instance that allows querying the graph and
other information (papers in a vector database, information from web APIs such
as OncoKB) via natural language. We provide two applications, one based on
[BioChatter Light](https://github.com/biocypher/biochatter-light), the other on
[BioChatter Next](https://github.com/biocypher/biochatter-next).

You can find more detailed information in a
[vignette](https://biochatter.org/vignettes/custom-decider-use-case/) on our
website and see hosted demonstrations of both applications at the following
links:

- Light: [https://decider-light.biocypher.org](https://decider-light.biocypher.org)

- Next: [https://decider-next.biocypher.org](https://decider-next.biocypher.org)

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
interface (for more info see this
[vignette](https://biochatter.org/vignettes/custom-bclight-simple/)). The Neo4j
instance can be accessed at `localhost:7474` and the BioChatter Light instance
at `localhost:8501`.

For the BioChatter Next variant, you can run the following command:

```{bash}
docker compose -f docker-compose-next.yml up -d
```

This will similarly build a KG, but also a Milvus instance for the vector
database, and instances of the `biochatter-server` REST API service and the
BioChatter Next app, which can be accessed at `localhost:3000`.

> [!IMPORTANT]
> For using OpenAI GPT as the language model, you will have to provide your API key as an environment variable (`OPENAI_API_KEY`) in your environment. You can do this using an export command (`export OPENAI_API_KEY=sk-...`) or by adding it to your bash profile; you could also provide it to Docker using an env file. We use GPT-3.5-turbo as the default model.

## Questions

The knowledge graph contains information about patients, genes, variants, drugs,
pathways, and clinical data. The schema of the KG can be seen below the query
interface as a JSON object. You can ask questions in natural language, such as:

- How many patients do we have, and what are their names?

- What was patient1's response to previous treatment, and which treatment did
they receive?

- Which patients have hr deficiency but have not received parp inhibitors?

- Does patient1 have a sequence variant in a gene that is druggable with
evidence level "1"? Which drug? Return unique values.

- Is there a patient with overlapping variants compared to patient1?

These are only few of infinitely many possible questions, and some may not
result in a valid query. The BioChatter Light interface allows manual
modification and rerunning of the query for prototyping and debugging.

## ⚙️ Local Installation

You can run the KG build locally using a virtual environment.

```{bash}
git clone https://github.com/biocypher/decider-genetics.git
cd decider-genetics
poetry install
poetry run python create_knowledge_graph.py
```
