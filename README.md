# COVID-Modeling
A simple SEIR for on/off-campus transmission of SARS-CoV-2 with delays in testing, behavioral modifiers, and infection/transmission modifiers. This stems from orginal work by Dr. Matt Ferrari at PSU . 

This uses very basic inputs to determine how many students are on-campus versus off-campus. We made some assumptions about student behavior on campus and added a "compliance" and "care.seeking" to parametrize how students will comply with isolation/contact-tracing and what percent of symptomatic students are expected to seek primary care. A certain number of tests, set by the user, are first allocated to classifying symptomatic tests. Extra tests are then disributed randomly to the rest of the student body.

This is NOT an individual-based model, so there are some odd work-arounds. The number of close contacts are randomly distributed with a poisson distribution with a mean/var equal to 7 here, though the end-user can change this for both on and off campus students. 

Testing delays are also incorportated here. Symptomatic cases are removed at the time of testing. Contacts are tested at the time that a close contact returns a positive test result, but removed after a delay. Asymptomatic individuals are tested randomly, and are removed only after a positive test result. So, with a 2 day delay in testing: a symptomatic person tested on day 0 would be removed on day 0. The close contacts of a symptomatic person would be contacted/tested on day 2. If positive, the close contacts are then removed on day 4. Asymptomatic individuals tested on day 0 would be removed after a positive test on day 2, etc. 

The vizualization file should be able to run on its own by calling the contents of the functions file. Let me know if you have issues there. 

This model is not peer-reviewed, nor do I try to claim any authority in epidemic modeling. I'm just an undergrad with a computer. Please use with skepticism and reach out if you have comments, ideas, or questions. 
