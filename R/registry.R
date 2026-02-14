epi_registry <- registry()

epi_registry$set_field("name", type = "character", is_key = TRUE)
epi_registry$set_field("model", type = "epi_model")
epi_registry$set_field("family", type = "character")
epi_registry$set_field("origin", type = "character")  # builtin | user
