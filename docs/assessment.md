# Assessment

How you are graded, how to submit, and what happens when things go wrong.

---

## Breakdown

| Component | Weight | Timing | Format |
|---|:---:|---|---|
| **Assignments (7, one per chapter)** | **20%** | Weeks 4, 6, 7, 10, 12, 14, 16 | Take-home; theory + code |
| **Mid-semester examination** | **30%** | Week 8 | 2 h, closed book |
| **Final examination** | **50%** | Week 18 | 3 h, closed book |
| **Total** | **100%** | | |

Per-assignment weight: **~3% each** (7 × 3% = 21%; the actual split is announced
with each assignment to reflect the size of the task — the total is exactly 20%).

---

## Assignments (20%)

### What to expect

Each assignment covers **one chapter** and combines theory with computation:

- **~60% theory.** Analytical problems similar to the tutorial problems, but
  longer and requiring a full solution — derivation, justification, and
  interpretation of the result.
- **~40% implementation.** A Python / Java / C# program that implements an
  algorithm from the chapter, runs it on a specified problem, and reports the
  results. Code must be **your own**; libraries such as `scipy.optimize` are
  allowed only when the assignment explicitly permits them.

### Rubric (per assignment)

| Criterion | Weight | What earns full marks |
|---|:---:|---|
| **Correctness of final answers** | 40% | All answers correct and complete, including units and interpretation |
| **Rigour of derivations** | 20% | Each step justified; no hidden assumptions; notation consistent |
| **Computational implementation** | 20% | Code runs, produces correct results, is readable and commented |
| **Interpretation and discussion** | 15% | Results placed in context; limitations and edge cases discussed |
| **Presentation** | 5% | Clear structure; plots labelled; references cited |

### Submission

- **Format.** A single `.pdf` report (typeset with LaTeX, Markdown → PDF, or Word)
  plus the source code in a `.zip` bundle. The PDF must **also** contain the
  code listings (or an appendix with the source).
- **Naming.** `A<n>-<lastname>-<student-id>.pdf` and `A<n>-<lastname>-<student-id>.zip`.
- **Channel.** Upload to the course Moodle/LMS assignment dropbox.
- **Deadline.** 23:59 local time on the due date. **Late submission loses 10%
  of the maximum mark per calendar day; after 3 days, the submission scores zero.**

---

## Mid-semester examination (30%)

- **Format.** 2 h, closed book, no calculators beyond a basic scientific one.
- **Content.** Covers **Chapters 1–3** in full. Weighting is approximately
  30% / 40% / 30% across the three chapters — reflecting the amount of material
  in each.
- **Structure.** Typically 4–6 questions, mixed theory and computation:
  state a definition or theorem; derive a result; solve a small problem by hand;
  interpret a numerical result.
- **Sample papers.** Past mid-semester papers are available on the LMS; the
  solutions are discussed in the tutorial in **week 7**.

### Common pitfalls

- **Confusing necessary and sufficient** conditions (FONC vs. SOSC).
- **Forgetting to check convexity** before claiming a stationary point is a global minimum.
- **Sign errors in the KKT multipliers** — the sign convention matters.
- **Wrong condition-number intuition** — a small condition number is *good*.
- **Copying a formula without stating the assumptions** — most formulas in
  Chapter 3 (penalty, barrier, augmented Lagrangian) require specific
  conditions (feasibility, strict feasibility, rank conditions).

---

## Final examination (50%)

- **Format.** 3 h, closed book, same calculator policy as the mid-semester.
- **Content.** Covers the **whole syllabus**, with approximately 40% on
  Chapters 1–3 and 60% on Chapters 4–7 (matching the split of teaching hours).
- **Structure.** Typically 6–8 questions. Theory and computation as in the
  mid-semester, with longer multi-part problems.
- **Emphasis.** Chapters 4–7 are examined more heavily because they contain most
  of the applied material — random walks, shot noise, cyclostationarity,
  sampling, spectral factorization, LMMSE, Wiener, Kalman, Markov chains, M/M/1,
  M/M/s.

### What to bring

- Non-programmable calculator.
- **One hand-written A4 formula sheet**, both sides. This must be *your own
  handwriting* — printed or photocopied sheets are not allowed.
- Student ID card.

The formula-sheet rule is deliberate: writing it forces you to review every
formula, and consulting it during the exam saves you the risk of a sign error.

---

## Grade boundaries

| Grade | Range | Description |
|---|---|---|
| **A** | 90–100 | Excellent — full understanding, no significant errors |
| **A−** | 85–89 | Very good — minor errors only |
| **B+** | 80–84 | Good — occasional errors, strong overall |
| **B** | 75–79 | Solid — some conceptual gaps but correct approach |
| **B−** | 70–74 | Adequate — noticeable gaps in understanding |
| **C+** | 65–69 | Fair — several conceptual errors |
| **C** | 60–64 | Satisfactory — meets minimum requirements |
| **F** | < 60 | Fail — does not meet requirements |

Final grades are **not curved**. The thresholds above are the standard ones for
the MSc program and are not adjusted after the fact.

---

## Late work and extensions

### Late policy

- **Assignments:** 10% of the maximum mark per calendar day (including weekends),
  capped at 3 days. After 3 days, the submission scores zero.
- **Examinations:** no late sitting. Missing an examination without documented
  cause scores zero.

### Extensions

Extensions are granted **only** for documented medical, family, or other
exceptional circumstances:

1. Submit a request **before** the deadline, with documentation
   (medical certificate, family emergency record).
2. The instructor will grant an extension of typically 3–7 days, proportional
   to the disruption.
3. If the deadline is missed due to a sudden emergency, contact the instructor
   within **48 hours** of the deadline with documentation.

Unsubstantiated requests are declined. Extensions are never granted retroactively
except for a confirmed emergency.

---

## Make-up examinations

- **Mid-semester:** students with documented medical or family emergencies
  may sit a **make-up examination in week 10**. The content and format are the
  same as the original; the maximum mark is capped at the class median of the
  original exam (standard practice to remove the "second chance" advantage).
- **Final:** governed by the university's rules for deferred examinations.
  Application is made to the Registrar's Office, with medical documentation,
  within 5 working days of the missed exam.

---

## Academic integrity

All submitted work (assignments, exams) must be **your own**.

### What is allowed

- Discussing the **concepts** and **methods** with classmates — this is
  encouraged.
- Consulting textbooks, lecture notes, and online references (cite them).
- Using standard libraries (`numpy`, `scipy`, `matplotlib`) as permitted by
  the assignment brief.
- Comparing your **final answers** with classmates **after** both have
  submitted — to check your reasoning.

### What is not allowed

- Copying code or text from another student, past submissions, or an online
  source.
- Sharing your code or write-up with classmates **before** the deadline.
- Using AI tools to generate code or text that you present as your own
  (see below).
- Re-submitting work submitted in another course.

### On AI tools

Large language models (ChatGPT, Claude, GitHub Copilot, etc.) are **powerful
learning tools** but **cannot be used to generate submitted work**. Permitted uses:

- Explaining a concept you don't understand.
- Debugging your own code (show the code, ask "why is this wrong?").
- Checking a derivation you have already done.

Not permitted: pasting an AI-generated answer or code as your submission.
If in doubt, ask the instructor **before** submitting.

### Consequences

Confirmed plagiarism or unauthorised AI use results in:

- **First offence:** zero on the assessment component, and a written warning on file.
- **Second offence:** zero on the assessment component and a report to the
  Dean's office, which may lead to academic probation or expulsion.

All submissions are screened with plagiarism-detection tools. Code is
additionally screened for AI-generated signatures.

---

## Feedback and grades

- **Assignment feedback:** returned within **2 weeks** of the submission
  deadline. Feedback is delivered through the LMS, with comments on each
  rubric criterion and a numeric mark.
- **Examination grades:** posted within **3 weeks** of the exam date.
- **Grade queries:** if you believe a mark is incorrect, email the instructor
  within **1 week** of receiving the result, with a written justification
  (what you expected, what you received, and why). Regrade requests that
  result in a lower mark are not reversed.

---

## Contact

- **Instructor:** Dr. Tesfamichael Gebregziabher
- **Email:** *(insert)*
- **Office:** *(insert)*
- **Office hours:** *(insert; e.g., Tue 14:00–16:00, Thu 10:00–12:00, or by appointment)*

For questions that would benefit the whole class (clarification of a syllabus
point, a typical exam question), post on the course forum — the instructor
monitors it daily and the answer helps everyone.
