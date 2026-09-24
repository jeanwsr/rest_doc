# ADR: Rules and Explanation

- **Status**: Proposed
- **Date**: 2026-09-16
- **Scope**: Authoring and maintenance of all ADRs in `rest_doc`; all developers participating in REST development.
- **Author**: GLM-5.3 (drafted in a mattpocock-skills:grill-with-docs session with Zhu Zhenyu (ajz34@outlook.com))
- **Reviewed by**: Zhu Zhenyu (ajz34@outlook.com)

## Decision

- REST records important technical decisions as ADRs: decisions that affect multiple modules or are hard to change later, together with their context, options, and rationale.
- Each ADR records exactly one decision, in a file named `adr-NNNN-slug.md`, numbered incrementally from ADR-0001; this document is the usage rules and explanation itself, not a numbered ADR.
- An ADR is a historical record, not a living document: once in effect, its substantive content is not rewritten; when a decision changes, a new ADR is written to take its place, and the old ADR receives only non-substantive updates — its status, or Supersedes/Partially superseded fields plus notices at the affected sections pointing to the successor (see Section 3.1) — leaving all other content unchanged.
- ADR status values: `Proposed` / `Rejected` / `In use` / `Disputed` / `Superseded (ADR-NNNN)` / `Deprecated`. A `Proposed` ADR is set to `In use` after maintainer confirmation; while `Disputed`, the decision remains binding until the dispute is resolved and otherwise disposed of.
- The Decision section of an ADR binds subsequent code within its scope; legacy code is allowed to deviate, which is stated in one sentence in the Scope field — no file-level conformance lists are maintained.
- This document is a living document: the rules may be revised directly by the maintainer with the date updated, without a successor ADR. The modification history is traceable via git blame.
- Writing or reviewing ADRs is led by human developers, optionally assisted by code agent skills such as `grill-with-docs` (see Section 4).

## 1. What is an ADR

An ADR (Architecture Decision Record) is a lightweight software-engineering practice: a short document recording one technical decision — what problem was faced at the time, which realistic alternatives existed, which one was chosen, and why.

"Architecture" does not mean that only grand system designs deserve an ADR. Anything that constrains future code over the long term is an architectural decision; this includes coding conventions, engineering techniques, testing processes, and so on.

The most important value of an ADR is recording the **why**. Code and API documentation can tell readers what a program *is* and how to use it, but cannot answer "why not the other approach that looks more natural". These reasons are known to everyone while the discussion is happening, yet a year later only guesses remain. By putting the reasoning, together with the rejected alternatives, on paper, future readers (including AI agents) can judge whether the premises of the time still hold today.

## 2. When to write an ADR: differences from other documents

| Document type | Question answered | Lifecycle | Examples |
| --- | --- | --- | --- |
| ADR | Why this decision was made | Written when the decision is made, kept as history | The documents in this directory |
| API documentation | How to use this interface/parameter | Continuously updated with the interface | rustdoc |
| Developer documentation | How the code is organized, how to contribute | Continuously maintained | Build guides, module descriptions |
| Enhancement proposal | What we intend to do in the future | Proposal—review—implementation, archived after implementation | Python PEP, NumPy NEP, Rust RFC |

The key difference is the lifecycle. API and developer documentation describe the **present** and must track the code continuously; enhancement proposals face the **future** and build consensus before implementation starts; an ADR is a **snapshot** recording the judgment at the moment the decision was made, and is not rewritten as the code evolves. The same piece of work often touches several document types at once: a new feature may start as an enhancement proposal for discussion, and after implementation its usage enters the API/user documentation; a technical choice, whether made before writing code or settled into a long-term convention during implementation, becomes an ADR.

REST currently has no formal proposal process or lifecycle maintenance: directional consensus is usually reached in person within the group or in Gitee issues/PRs, and often simply lands as an implementation. The process described above is therefore common industry practice rather than REST's current state.

Consider writing an ADR when any of the following holds:

- The decision's constraints extend beyond, or will extend beyond, a single module or a single PR;
- The decision is hard to reverse, or expensive to reverse, once made;
- The decision is a trade-off among several realistic options, and it is foreseeable that someone will ask "why not that way".

Conversely, local implementation details, bug fixes, and refactors that introduce no new decision do not need an ADR; code comments or PR descriptions suffice.

## 3. How to write an ADR

### 3.1 Structure

ADRs live in the `contributor/adr/` directory of each language tree (`rest_doc/source_zh/contributor/adr/` for Chinese, `rest_doc/source/contributor/adr/` for English), in files named `adr-NNNN-slug.md`, and are registered in the toctree of that directory's `index.md`. The structure is:

- Header metadata, one item per line:
  - **Status**: `Proposed` / `Rejected` / `In use` / `Disputed` / `Superseded (ADR-NNNN)` / `Deprecated`;
  - **Date**: the date the decision actually took effect. If the decision landed together with a code mechanism, take the date the mechanism landed (verifiable from git history), not the day the document was written;
  - **Scope**: which crates/modules the decision constrains; deviations in legacy code are stated here in one sentence;
  - **Author**: the drafter, either a human developer or an AI agent;
  - **Reviewed by** (optional): for agent-drafted records, the human developer responsible for confirmation;
  - **Supersedes** (optional): the older ADR this record takes over, either the whole record or a specific section (e.g. `ADR-0001 Section 3.3`);
  - **Partially superseded** (optional): for an in-use ADR of which only individual parts are taken over by a newer ADR, list the superseded sections and the successor; a dated one-line notice pointing to the successor may be placed under the affected section heading, while the original body content is kept unchanged.
- **Decision**: the normative content. Compressed into a bulleted list of a few immediately readable points, optionally with one minimal code example; most readers most of the time read only this part, so it must stand alone, independent of the body.
- **Body**: narrative sections — design motivation, usage, technical details, alternatives, remarks, etc.; the body is explanatory (non-normative) by default, and individual sections may be marked otherwise (e.g. "non-normative"). The body is what distinguishes an ADR from a dry rule list: only when the reasoning is spelled out can readers judge for themselves whether the premises have changed.

A status of `Deprecated` means the whole decision is withdrawn with no successor, which is a different thing from in-text markers on individual rejected alternatives (see Section 3.2).

### 3.2 Writing conventions

- Rejected or deprecated alternatives are **kept in the text**, clearly marked with wording such as "deprecated"; they are the best evidence of "why not".
- No file-level conformance lists (which files comply, which do not). Such lists rot quickly as the code evolves; a single sentence "legacy code partially does not comply yet and will be aligned gradually during refactoring" suffices.
- When pointing to code or external material, give stable locations (module paths, source links); avoid phrasings that rot over time, such as "the current PR".
- For writing or translation into Simplified Chinese, code, function names, and proper nouns are kept in the original language.

## 4. Workflow in the AI era

Most REST developers are not professional software engineers; until one can build ADRs independently, it is recommended to work with code agent skills such as `grill-with-docs` (from [mattpocock/skills](https://github.com/mattpocock/skills)) or similar tools. Two working modes are recommended:

- **Developer-led prompting, AI-led ADR drafting**: the developer and the agent go through rounds of Q&A. The agent reads the relevant code, verifies facts, and asks sharp questions about each branch of the design; the decision authority always stays with the developer. Once the design is clarified, the agent drafts the ADR and the developer edits and finalizes it.
- **Developer writes the ADR, AI reviews**: for a developer-written ADR, the agent performs fact analysis and structural adjustment. Fact analysis means the agent checks every claim against the current code and points out outdated, ambiguous, or insufficiently justified statements. Structural adjustment means the agent reorganizes the content so that it better conforms to ADR norms and readability, and better supports collaboration with AI agents.

Emphasizing human-developer leadership:
- ADR documents are strongly binding. A wrong ADR causes more code drift than having no ADR at all.
- AI agents can produce fluent text quickly, but unverified claims within it read just as fluently; the round-by-round Q&A process forces every claim written into an ADR through fact-checking and developer confirmation.

ADRs are also documents written for AI agents: an agent that reads the ADRs before generating or modifying code can follow their norms and taboos instead of stepping into pits that were already stepped into.

**It is not recommended to draft the ADR and implement the code within the same AI agent session.** The output of an AI agent is strongly influenced by the prompt. Doing both in one session tends to make the agent, while writing the ADR, steer toward implementation details suggested by the prompt, drifting away from the intent of the ADR.

## 5. Further reading

- [architecture-decision-record](https://github.com/architecture-decision-record/architecture-decision-record): a brief introduction to ADRs and consolidated resources (templates, tools, examples).
