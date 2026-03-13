# Reference Outputs

`nlswork_stata.toml` is the checked-in Stata baseline used by the maintenance workflow.

If you want to compare against your boss''s output as well, add a new file named `boss_reference.toml` in this folder using the same schema. The maintenance tests automatically load every `*.toml` file here except `*.example.toml`.
