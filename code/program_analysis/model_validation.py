def validate_model_object(model):
    print(type(model))

    print("\nNetwork:")
    print(type(model.network))
    print(model.network.shape)

    print("\nResponses:")
    print(type(model.responses))
    print(model.responses.shape)

    print("\nRelative responses:")
    print(type(model.relative_responses))
    print(model.relative_responses.shape)

    print("\nSample list:")
    print(model.sample_list)    
    return 