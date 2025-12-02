# bezier_utils/utils/event_utils.py

def updateMetaBtns(caller, event, keymap = None):
    if(keymap is None):
        keymap = {'LEFT_SHIFT': 'shift', 'RIGHT_SHIFT': 'shift',
            'LEFT_CTRL':'ctrl', 'RIGHT_CTRL':'ctrl',
            'LEFT_ALT': 'alt', 'RIGHT_ALT': 'alt'}

    var = keymap.get(event.type)

    if(var is not None):
        expr = 'caller.' + var + ' = '
        if(event.value == 'PRESS'): exec(expr +'True')
        if(event.value == 'RELEASE'): exec(expr +'False')
        return True

    return False
