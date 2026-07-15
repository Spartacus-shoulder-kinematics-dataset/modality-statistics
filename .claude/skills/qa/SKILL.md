---
name: qa
description: Document the recent conversation into Q&A.md as a synthetic question/answer digest so the user can revisit and digest the reasoning. Invoke with /qa (the user may call it "/Q&A"). Optionally pass a topic to focus on.
---

# Q&A knowledge digest

Goal: help the user **digest the knowledge** from our interactions by turning the
recent conversation into a compact question/answer log in `Q&A.md` at the project
root.

When invoked:

1. **Read `Q&A.md`** if it exists (create it with the header below if not). Find
   the last logged entry so you only add what is new. If the user passed a topic
   argument, focus the entry on that topic.

2. **Reformulate the user's question(s)** since the last entry into one concise,
   self-contained question each — strip the back-and-forth, keep the intent.

3. **Answer synthetically** — 2–5 lines per question. State what was decided or
   done and *why*, in digest form (not a transcript). Link to the concrete
   artifact when useful (`analysis_shoulder.R`, `figures/...`, a file path).

4. **At every decision point, propose exactly THREE options.** Mark the chosen
   one with ★ and put it first; give the two alternatives with a one-line
   trade-off each. This is mandatory whenever a choice was (or should be) made —
   it is the core of this skill.

5. **Append** the new entry/entries to `Q&A.md` in chronological order, numbered,
   each dated (use the current date). Keep it terse — this is a study aid, not a
   report. Do not rewrite past entries.

## Entry format

```markdown
## NN. <Reformulated question>   _(YYYY-MM-DD)_

**A.** <synthetic answer, 2–5 lines, with *why*.>

**Decision — <what was decided>:**
- ★ **<chosen option>** — <one-line reason it won>
- <alternative 2> — <one-line trade-off>
- <alternative 3> — <one-line trade-off>
```

Omit the Decision block only when the question genuinely involved no choice.
Include one Decision block per real decision (an entry may have several).

## File header (create once, at the top of Q&A.md)

```markdown
# Q&A — decisions & reasoning digest

Synthetic log of what was asked, what was decided, and the alternatives — kept so
the reasoning can be revisited. Generated/updated by the `/qa` skill.
```
